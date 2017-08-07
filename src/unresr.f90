module unresm
   ! module to provide unresr for NJOY2016
   use locale
   implicit none
   private
   public unresr

   !global variables for unresr
   integer::nendf,nin,nout,nscr
   integer,parameter::nmtx=5
   integer::nunr,matd,nsigz,mtx(nmtx),npnts(nmtx)
   integer::intunr,lssf,nro,naps
   integer::isot(20),modet(20),ibaset(20),nsect
   real(kr)::abnt(20),elt(20),eht(20)
   real(kr)::tr(62,62),ti(62,62)
   real(kr),dimension(:),allocatable::tempu

contains

   subroutine unresr
   !--------------------------------------------------------------------
   !
   ! compute unresolved resonance cross-sections
   !
   ! The method of ETOX is used to compute self-shielded
   ! unresolved resonance cross-sections on the energy grid of
   ! the unresolved parameters.  Subsequent interpolation is
   ! to be on the cross-sections and not on the parameters.
   ! Additional energy grid points are added at quarter lethargy
   ! intervals if only three or fewer grid points are found.
   ! The accurate Hwang quadrature set is used for the integrals.
   !
   !---input specifications (free format)---------------------------
   !
   ! card 1
   !   nendf   unit for endf tape
   !   nin     unit for input pendf tape
   !   nout    unit for output pendf tape
   ! card 2
   !   matd    material to be processed
   !   ntemp   no. of temperatures (default=1)
   !   nsigz   no. of sigma zeroes (default=1)
   !   iprint  print option (0=min, 1=max) (default=0)
   ! card 3
   !   temp    temperatures in Kelvin (including zero)
   ! card 4
   !   sigz    sigma zero values (including infinity)
   !       cards 2, 3, 4 must be input for each material desired
   !       matd=0/ terminates execution of unresr.
   !
   !--------------------------------------------------------------------
   use mainio ! provides nsysi,nsyso
   use endf   ! provides endf routines and variables
   use util   ! provides timer,error,openz,repoz,closz
   ! locals
   integer::i,it,nb,nw,iprint,ntemp,lrp,l,ie
   integer::ix,is,j,new,nx,nc,nd,newmat,ncds,mfd,mtd
   integer::k,nwd,lstart,nwds
   real(kr)::time,temz,za,awr,ez
   real(kr)::bkgz(4)
   character(60)::strng1,strng2
   real(kr),dimension(:),allocatable::scr
   real(kr),dimension(:),allocatable::c,d
   real(kr),dimension(:),allocatable::eunr
   real(kr),dimension(:,:),allocatable::sb
   real(kr),dimension(:),allocatable::arry
   real(kr),dimension(:),allocatable::b
   real(kr),dimension(:),allocatable::temp
   real(kr),dimension(:),allocatable::sigz
   real(kr),dimension(:,:),allocatable::sigu
   integer,parameter::maxscr=1000
   integer,parameter::maxeunr=150
   integer,parameter::maxarry=10000

   !--initialize
   nscr=10
   nsh=0
   mtx(1)=1
   mtx(2)=2
   mtx(3)=18
   mtx(4)=19
   mtx(5)=102
   call timer(time)
   read(nsysi,*) nendf,nin,nout
   if (nin.lt.0.and.nout.gt.0) call error('unresr',&
      'mode conversion between nin and nout not allowed.',' ')
   if (nin.gt.0.and.nout.lt.0) call error('unresr',&
      'mode conversion between nin and nout not allowed.',' ')
   if (nin.lt.0) nscr=-nscr
   call openz(nendf,0)
   call openz(nin,0)
   call openz(nout,1)
   call openz(nscr,1)
   write(nsyso,'(/'' unresr...'',&
     &''calculation of unresolved resonance cross sections'',&
     &9x,f8.1,''s'')') time
   write(nsyse,'(/'' unresr...'',59x,f8.1,''s'')')  time
   allocate(scr(maxscr))
   write(nsyso,'(/&
     &'' unit for input endf tape ............. '',i10/&
     &'' unit for input pendf tape ............ '',i10/&
     &'' unit for output pendf tape ........... '',i10)')&
     nendf,nin,nout
   call uwtab
   call repoz(nin)
   call repoz(nendf)
   call repoz(nout)
   call tpidio(nin,nout,0,scr,nb,nw)
   call tpidio(nendf,0,0,scr,nb,nw)

   !--loop over requested materials
  110 continue
   iprint=0
   ntemp=1
   nsigz=1
   read(nsysi,*) matd,ntemp,nsigz,iprint
   if (matd.eq.0) go to 400
   if (allocated(temp)) then
      deallocate(temp)
      deallocate(sigz)
      deallocate(sigu)
   endif
   allocate(temp(ntemp))
   allocate(sigz(nsigz))
   allocate(sigu(5,nsigz))
   read(nsysi,*) (temp(i),i=1,ntemp)
   read(nsysi,*) (sigz(i),i=1,nsigz)
   write(nsyso,'(/&
     &'' temperatures ......................... '',1p,e10.3/&
     &(40x,e10.3))') (temp(i),i=1,ntemp)
   write(nsyso,'(&
     &'' sigma zero values .................... '',1p,e10.3/&
     &(40x,e10.3))') (sigz(i),i=1,nsigz)
   write(nsyso,'(&
     &'' print option (0 min., 1 max.) ........ '',i10)') iprint

   !--loop over requested temperatures
   do 320 it=1,ntemp
   temz=temp(it)
   call timer(time)
   write(nsyso,'(/&
     &'' mat = '',i4,3x,'' temp = '',1p,e10.3,37x,0p,f8.1,''s'')')&
     matd,temz,time

   !--read resonance parameters from file 2 on endf tape.
   !--read background cross sections from file 3
   do i=1,nmtx
      npnts(i)=0
   enddo
   if (it.gt.1) go to 130
   nunr=0
   call findf(matd,1,451,nendf)
   call contio(nendf,0,0,scr,nb,nw)
   lrp=l1h
   if (lrp.eq.0) go to 340
   call contio(nendf,0,0,scr,nb,nw)
   if (n1h.ne.0) then
      iverf=4
   else if (n2h.eq.0) then
      iverf=5
   else
      iverf=6
   endif
   call findf(matd,2,151,nendf)
   call contio(nendf,0,0,scr,nb,nw)
   za=c1h
   awr=c2h
   allocate(eunr(maxeunr))
   allocate(arry(maxarry))
   call rdunf2(eunr,scr,arry,maxarry)
   if (nunr.eq.0) go to 350
   if (lssf.gt.0) go to 130
   call findf(matd,3,1,nendf)
   allocate(sb(nunr,4))
   call rdunf3(eunr,sb,nunr,scr)
  130 continue
   nb=12+nsigz+nunr*(1+5*nsigz)
   allocate(b(nb))
   l=1
   b(l)=za
   b(l+1)=awr
   b(l+2)=lssf
   b(l+3)=0
   b(l+4)=0
   b(l+5)=intunr
   b(l+6)=temz
   b(l+7)=0
   b(l+8)=5
   b(l+9)=nsigz
   b(l+10)=nsigz+nunr*(1+5*nsigz)
   b(l+11)=nunr
   l=l+11
   ! store sigz into next nsigz location of a.
   do i=1,nsigz
      l=l+1
      b(l)=sigz(i)
   enddo

   !--compute unresolved resonance cross-sections
   !--for all grid points and values of sigma zero.
   do ie=1,nunr
      ez=eunr(ie)
      l=l+1
      b(l)=ez
      if (lssf.le.0) then
         bkgz(1)=sb(ie,1)
         bkgz(2)=sb(ie,2)
         bkgz(3)=sb(ie,3)
         bkgz(4)=sb(ie,4)
      else
         bkgz(1)=0
         bkgz(2)=0
         bkgz(3)=0
         bkgz(4)=0
      endif
      if (iprint.eq.1) write(nsyso,&
        '('' energy ='',1p,e12.4)') ez
      call unresl(ez,temz,sigz,nsigz,bkgz,sigu,arry)
      do ix=1,5
         do is=1,nsigz
            l=l+1
            b(l)=sigfig(sigu(ix,is),7,0)
         enddo
         if (l.gt.nb)&
           call error('unresr','storage exceeded.',' ')
         if (iprint.eq.1) write(nsyso,'(1x,1p,8e11.3)')&
           (sigu(ix,j),j=1,nsigz)
      enddo
   enddo
   if (nin.eq.0) go to 310
   if (nout.eq.0) go to 310

   !--check for previous mt152
   new=1
   call findf(matd,2,151,nin)
   call tosend(nin,0,0,scr)
   call contio(nin,0,0,scr,nb,nw)
   if (mth.eq.152) new=0
   call findf(matd,1,451,nin)

   !--write new pendf tape
   call contio(nin,0,0,scr,nb,nw)
   if (iverf.lt.5) then
      nx=n2h
      if (nx.gt.0.and.new.gt.0) scr(6)=scr(6)+1
   endif
   call contio(0,nout,0,scr,nb,nw)
   if (iverf.ge.5) call contio(nin,nout,0,scr,nb,nw)
   if (iverf.eq.6) call contio(nin,nout,0,scr,nb,nw)
   call hdatio(nin,0,0,scr,nb,nw)
   scr(1)=temz
   if (iverf.ge.5) then
      nx=n2h
      if (nx.gt.0.and.new.gt.0) scr(6)=scr(6)+1
   endif
   call hdatio(0,nout,0,scr,nb,nw)
   do while (nb.ne.0)
      call moreio(nin,nout,0,scr,nb,nw)
   enddo
   nw=nx
   if (nw.ne.0) then
      nc=6*nx+6
      allocate(c(nc))
      nd=6*nx
      allocate(d(nd))
      call dictio(nin,0,0,d,nb,nw)
      j=0
      newmat=0
      ncds=2+(nint(b(11))-1)/6
      ix=1
      do i=1,nw,6
         mfd=nint(d(i+2))
         mtd=nint(d(i+3))
         if (mfd.gt.2.or.mtd.eq.152) then
            if (newmat.le.0) then
               newmat=1
               c(j+1)=0
               c(j+2)=0
               c(j+3)=2
               c(j+4)=152
               c(j+5)=ncds
               c(j+6)=0
               if (new.gt.0) j=j+6
            endif
         endif
         do k=1,6
            c(k+j)=d(k+i-1)
         enddo
         j=j+6
      enddo
      nwd=j/6
      call dictio(0,nout,0,c,nb,nwd)
      deallocate(d)
      deallocate(c)
   endif
   ! copy through file 2, mt 151
   call tofend(nin,nout,0,scr)
   call tosend(nin,nout,0,scr)
   if (new.eq.0) call tosend(nin,0,0,scr)
   ! write new mt on output tape
   math=matd
   mfh=2
   mth=152
   call contio(0,nout,0,b,nb,nw)
   lstart=7
   nwds=l-6
   nw=nwds
   if (nw.gt.npage+6) nw=npage+6
   call listio(0,nout,0,b(lstart),nb,nw)
  210 continue
   lstart=lstart+nw
   nwds=nwds-nw
   nw=nwds
   if (nw.gt.npage) nw=npage
   if (nb.eq.0) go to 220
   call moreio(0,nout,0,b(lstart),nb,nw)
   go to 210
  220 continue
   call asend(nout,0)
   ! copy file end card and rest of material
   call tomend(nin,nout,0,scr)

   !--write report of calculation
  310 continue
   call timer(time)
   write(nsyso,'('' generated cross sections at '',i3,&
     &'' points'',30x,f8.1,''s'')') nunr,time
   deallocate(b)
  320 continue
   if (lssf.le.0) deallocate(sb)
   deallocate(arry)
   deallocate(eunr)
   go to 110

   !--mat has no resonance parameters. copy as is to nout.
  340 continue
   write(strng1,&
     '(''mat'',i5,'' has no resonance parameters'')') matd
   write(strng2,'(''copy as is to nout'')')
   call mess('unresr',strng1,strng2)
   go to 355

   !--mat has no unresolved.  copy as is to nout
  350 continue
   write(strng1,&
     '(''mat'',i5,'' has no unresolved parameters'')') matd
   write(strng2,'(''copy as is to nout'')')
   call mess('unresr',strng1,strng2)
  355 continue
   call findf(matd,1,451,nin)
    it=0
  360 continue
   call contio(nin,0,0,scr,nb,nw)
   if (math.ne.matd) go to 110
   call contio(0,nout,0,scr,nb,nw)
   call tomend(nin,nout,0,scr)
   it=it+1
   if (it.lt.ntemp) go to 360
   if (allocated(eunr)) deallocate(eunr)
   if (allocated(arry)) deallocate(arry)
   go to 110

   !--finished with unresr.
  400 continue
   call atend(nout,0)
   call repoz(nout)
   call repoz(nin)
   call repoz(nendf)
   call closz(nendf)
   call closz(nin)
   call closz(nout)
   call closz(nscr)
   if (allocated(temp)) then
      deallocate(temp)
      deallocate(sigz)
      deallocate(sigu)
   endif
   if (allocated(tempu)) deallocate(tempu)
   call timer(time)
   write(nsyso,'(/69x,f8.1,''s''/&
     &1x,7(''**********''),''*******'')') time
   return
   end subroutine unresr

   subroutine rdunf2(eunr,scr,arry,jx)
   !--------------------------------------------------------------------
   ! Read and/or copy unresolved resonance parameters from file 2.
   !--------------------------------------------------------------------
   use endf ! provides endf routines and variables
   use util ! provides error,mess,sigfig
   ! externals
   integer::jx
   real(kr)::eunr(*)
   real(kr)::scr(*)
   real(kr)::arry(*)

   ! internals
   integer::indep,nis,inow,is,nb,nw,lfw,ner,ier,lru,lrf,nls,l
   integer::njs,j,mode,ll,jnow,n,i,ne,inow1,ir,k,ist,loc
   integer::ii,jj,iovlp,lim
   real(kr)::elr,ehr,abn,el,eh,et,spi,ay,awri,enow,elast,enext
   real(kr)::enut,en
   character(60)::strng
   integer,parameter::ngridu=78
   real(kr),dimension(78)::egridu=(/&
     1.0e1_kr,1.25e1_kr,1.5e1_kr,1.7e1_kr,2.0e1_kr,2.5e1_kr,3.0e1_kr,&
     3.5e1_kr,4.0e1_kr,5.0e1_kr,6.0e1_kr,7.2e1_kr,8.5e1_kr,1.0e2_kr,&
     1.25e2_kr,1.5e2_kr,1.7e2_kr,2.0e2_kr,2.5e2_kr,3.0e2_kr,3.5e2_kr,&
     4.0e2_kr,5.0e2_kr,6.0e2_kr,7.2e2_kr,8.5e2_kr,1.0e3_kr,1.25e3_kr,&
     1.5e3_kr,1.7e3_kr,2.0e3_kr,2.5e3_kr,3.0e3_kr,3.5e3_kr,4.0e3_kr,&
     5.0e3_kr,6.0e3_kr,7.2e3_kr,8.5e3_kr,1.0e4_kr,1.25e4_kr,1.5e4_kr,&
     1.7e4_kr,2.0e4_kr,2.5e4_kr,3.0e4_kr,3.5e4_kr,4.0e4_kr,5.0e4_kr,&
     6.0e4_kr,7.2e4_kr,8.5e4_kr,1.0e5_kr,1.25e5_kr,1.5e5_kr,1.7e5_kr,&
     2.0e5_kr,2.5e5_kr,3.0e5_kr,3.5e5_kr,4.0e5_kr,5.0e5_kr,6.0e5_kr,&
     7.2e5_kr,8.5e5_kr,1.e6_kr,1.25e6_kr,1.5e6_kr,1.7e6_kr,2.0e6_kr,&
     2.5e6_kr,3.0e6_kr,3.5e6_kr,4.0e6_kr,5.0e6_kr,6.0e6_kr,7.2e6_kr,&
     8.5e6_kr/)
   real(kr),parameter::wide=1.26e0_kr
   integer,parameter::n150=150
   real(kr),parameter::onemev=1.e6_kr
   real(kr),parameter::small=1.e-12_kr
   real(kr),parameter::step=1.01e0_kr
   real(kr),parameter::zero=0

   !--initialize
   indep=0
   elr=0
   ehr=onemev
   nunr=1
   eunr(1)=onemev

   !--head card read in calling program
   nis=n1h
   nsect=0
   inow=1

   !--loop over isotopes
   do 300 is=1,nis
   call contio(nendf,0,0,scr,nb,nw)
   abn=c2h
   lfw=l2h
   ner=n1h
   ier=0

   !--loop over all energy ranges
  110 continue
   call contio(nendf,0,0,scr,nb,nw)
   ier=ier+1
   el=sigfig(c1h,7,0)
   eh=sigfig(c2h,7,0)
   lru=l1h
   lrf=l2h
   nro=n1h
   naps=n2h
   if (lru.eq.2) go to 140
   if (eh.gt.elr) elr=eh

   !--read through this sub-section
   call contio(nendf,0,0,scr,nb,nw)
   nls=n1h
   if (lrf.eq.4.or.lrf.eq.7) then
      call listio(nendf,0,0,scr,nb,nw)
      do while (nb.ne.0)
         call moreio(nendf,0,0,scr,nb,nw)
      enddo
   endif
   do l=1,nls
      if (lrf.ne.4) then
         call listio(nendf,0,0,scr,nb,nw)
         do while (nb.ne.0)
            call moreio(nendf,0,0,scr,nb,nw)
         enddo
         if (lrf.eq.7) then
            call listio(nendf,0,0,scr,nb,nw)
            do while (nb.ne.0)
               call moreio(nendf,0,0,scr,nb,nw)
            enddo
         endif
      else
         call contio(nendf,0,0,scr,nb,nw)
            njs=n1h
         do j=1,njs
            call listio(nendf,0,0,scr,nb,nw)
            do while (nb.ne.0)
            call moreio(nendf,0,0,scr,nb,nw)
            enddo
         enddo
      endif
   enddo
   if (ier.lt.ner) go to 110
   go to 300

   !--read and copy unresolved resonance parameters
  140 continue
   nsect=nsect+1
   isot(nsect)=is
   abnt(nsect)=abn
   elt(nsect)=el
   eht(nsect)=eh
   mode=lrf+4*(lru-1)
   modet(nsect)=mode
   ibaset(nsect)=inow
   if (eh.lt.ehr) ehr=eh
   ! add el and eh to list of nodes
   ! shade nodes around discontinuities
   et=sigfig(el,7,-1)
   call ilist(et,eunr,nunr)
   et=sigfig(el,7,+1)
   call ilist(et,eunr,nunr)
   et=sigfig(eh,7,-1)
   call ilist(et,eunr,nunr)
   et=sigfig(eh,7,+1)
   call ilist(et,eunr,nunr)

   ! if present, read and store the energy-dependent
   ! scattering radius tab1 data.
   if (nro.eq.1) then
      call tab1io(nendf,0,0,arry(inow),nb,nw)
      jj=inow+nw
      do while (nb.ne.0)
         call moreio(nendf,0,0,arry(jj),nb,nw)
         jj=jj+nw
         if (jj.gt.jx) then
            call error('rdunf2','storage in a exceeded',' ')
         endif
      enddo
      inow=jj
   endif

   !--branch to specified representation
   if (lrf.eq.2) go to 250
   intunr=2
   if (lfw.eq.1) go to 190

   !--all parameters independent of energy
   call contio(nendf,0,0,scr,nb,nw)
   spi=c1h
   ay=c2h
   lssf=l1h
   nls=n1h
   ! loop over l states
   do l=1,nls
      call listio(nendf,0,0,scr,nb,nw)
      awri=c1h
      ll=l1h
      njs=n2h
      if (l.eq.1) then
         arry(inow)=awri
         arry(inow+1)=ay
         arry(inow+2)=spi
         arry(inow+3)=lfw
         arry(inow+4)=nls
         arry(inow+5)=0
         inow=inow+6
      endif
      arry(inow)=njs
      arry(inow+1)=ll
      inow=inow+2
      jnow=6
      do n=1,njs
         do i=1,5
            arry(inow)=scr(jnow+i)
            inow=inow+1
         enddo
         jnow=jnow+6
      enddo
   enddo
   if (inow.gt.jx) call error('rdunf2','storage exceeded.',' ')
   indep=1
   go to 290

   !--fission widths energy dependent
  190 continue
   call listio(nendf,0,0,scr,nb,nw)
   spi=c1h
   ay=c2h
   lssf=l1h
   ne=n1h
   nls=n2h
   arry(inow+1)=ay
   arry(inow+2)=spi
   arry(inow+3)=lfw
   arry(inow+4)=nls
   arry(inow+5)=ne
   ! store fission width energy points
   inow1=inow+5
   do i=1,ne
      enow=sigfig(scr(i+6),7,0)
      arry(inow1+i)=enow
      ! add to list of energy nodes
      if (i.gt.1.and.i.le.ne) then
         call ilist(enow,eunr,nunr)
      endif
   enddo
   ! loop over l states
   do l=1,nls
      call contio(nendf,0,0,scr,nb,nw)
      if (l.eq.1) then
         arry(inow)=c1h
         inow=inow+6+ne
      endif
      njs=n1h
      arry(inow)=njs
      arry(inow+1)=l1h
      inow=inow+2
      ! loop over j states
      do j=1,njs
         call listio(nendf,0,0,scr,nb,nw)
         ne=n1h-6
         arry(inow)=l2h
         ! store parameters
         do ir=1,5
            arry(inow+ir)=scr(6+ir)
         enddo
         inow1=inow+5
         ! store fission widths
         do k=1,ne
            arry(k+inow1)=scr(k+12)
         enddo
         inow=inow+6+ne
      enddo
   enddo
   if (inow.gt.jx) call error('rdunf2','storage exceeded.',' ')
   go to 290

   !--all parameters energy dependent
  250 continue
   ist=inow
   ! first cont record
   call contio(nendf,0,0,scr,nb,nw)
   lssf=l1h
   nls=n1h
   arry(inow+1)=c2h
   arry(inow+2)=c1h
   arry(inow+3)=n1h
   inow=inow+4
   ! loop over l states
   do l=1,nls
      call contio(nendf,0,0,scr,nb,nw)
      if (l.eq.1) arry(ist)=c1h
      njs=n1h
      arry(inow)=njs
      arry(inow+1)=l1h
      inow=inow+2
      ! loop over j states
      do j=1,njs
         call listio(nendf,0,0,scr,nb,nw)
         loc=1+nw
         do while (nb.ne.0)
            call moreio(nendf,0,0,scr(loc),nb,nw)
            loc=loc+nw
         enddo
         ne=n2h
         arry(inow)=c1h
         arry(inow+1)=l1h
         intunr=l1h
         arry(inow+2)=ne
         ! store numbers of degrees of freedom
         do k=3,6
            arry(k+inow)=scr(k+6)
         enddo
         inow=inow+7
         jnow=12
         ! store parameters
         do n=1,ne
            do k=1,6
               if (scr(jnow+k).le.zero) scr(jnow+k)=small
               arry(inow+k-1)=scr(jnow+k)
            enddo
            inow=inow+6
            ! add to list of energy nodes
            if (n.ne.1.and.n.ne.ne.and.l.eq.1.and.j.eq.1) then
               enow=sigfig(scr(jnow+1),7,0)
               call ilist(enow,eunr,nunr)
            endif
            jnow=jnow+6
         enddo
      enddo
   enddo
   if (inow.gt.jx) call error('rdunf2','storage exceeded.',' ')

   !--continue loop over isotopes and energy ranges
  290 continue
   if (ier.lt.ner) go to 110
  300 continue
   nunr=nunr-1
   if (nunr.eq.0) return

   !--add extra nodes if needed
   !--remove first and last energy nodes
   !--flag resolved-unresolved overlap energies
   nunr=nunr+1
   i=1
   elast=eunr(2)
  350 continue
   i=i+1
   enext=eunr(i+1)
   if (enext.ge.onemev) go to 380
   if (enext.lt.wide*elast.and.indep.eq.0) go to 370
   et=elast
  360 continue
   do 362 ii=1,ngridu
   enut=egridu(ii)
   if (enut.gt.step*et) go to 363
  362 continue
   enut=enext
  363 continue
   et=enut
   if (et.ge.enext) go to 370
   call ilist(et,eunr,nunr)
   i=i+1
   go to 360
  370 continue
   elast=eunr(i+1)
   go to 350
  380 continue
   iovlp=0
   lim=nunr-1
   nunr=0
   en=0
   do i=2,lim
      et=eunr(i)
      if (et.ge.en) then
         if (et.lt.elr) et=-et
         if (et.lt.zero) iovlp=1
         nunr=nunr+1
         eunr(nunr)=et
         en=sigfig(abs(et),7,2)
      endif
   enddo
   nunr=nunr-1
   if (iovlp.eq.1) call mess('rdunf2',&
     'resolved-unresolved overlap energies',&
     'marked with minus signs.')
   if (eunr(nunr).gt.ehr) then
      write(strng,&
        '(''unresolved-smooth overlap above e='',1p,e12.4)') ehr
      call mess('rdunf2',strng,' ')
   endif
   return
   end subroutine rdunf2

   subroutine ilist(e,elist,nlist)
   !--------------------------------------------------------------------
   ! Insert a new energy into a list of enegies in ascending order.
   ! Omit duplicate values.  Initial list must be primed with one
   ! energy larger than any others to be added.
   !--------------------------------------------------------------------
   ! externals
   integer::nlist
   real(kr)::e,elist(*)
   ! internals
   integer::jl,i,j

   jl=0
   do 110 i=1,nlist
   if (e.gt.elist(i)) go to 110
   if (e.eq.elist(i)) go to 130
   jl=nlist-i+1
   do j=1,jl
      elist(nlist-j+2)=elist(nlist-j+1)
   enddo
   elist(i)=e
   go to 130
  110 continue
   jl=1
   elist(nlist+1)=e
  130 continue
   if (jl.gt.0) nlist=nlist+1
   return
   end subroutine ilist

   subroutine rdunf3(eunr,sb,nunr,scr)
   !--------------------------------------------------------------------
   ! Read unresolved region background cross sections from file 3.
   ! Compute the background cross sections on the unresolved
   ! energy grid and store for use in unresl.
   !--------------------------------------------------------------------
   use endf ! provides contio,gety1

   ! external
   integer::nunr
   real(kr)::eunr(nunr),sb(nunr,4),scr(*)
   ! internals
   integer::nb,nw,ix,kx,ie,idis,np
   real(kr)::sig,e,enext,bkg
   real(kr),parameter::up=1.00001e0_kr
   real(kr),parameter::dn=0.99999e0_kr

   ix=1
  100 continue
   call contio(nendf,0,0,scr,nb,nw)
  105 continue
   kx=ix
   if (mtx(ix).gt.18) kx=ix-1
   if (mth.lt.mtx(ix)) go to 120
   if (mth.eq.mtx(ix)) go to 130
   do ie=1,nunr
      sb(ie,kx)=0
   enddo
   go to 160
  120 continue
   call tosend(nendf,0,0,scr)
   go to 100
  130 continue
   if (mtx(ix).eq.19) go to 150
   e=0
   call gety1(e,enext,idis,bkg,nendf,scr)
   np=nint(scr(6))
   npnts(ix)=np
   do ie=1,nunr
      e=abs(eunr(ie))
      if (ie.eq.1) e=up*e
      if (ie.eq.nunr) e=dn*e
      call gety1(e,enext,idis,sig,nendf,scr)
      sb(ie,kx)=sig
   enddo
   call tosend(nendf,0,0,scr)
   go to 160
  150 continue
   npnts(ix)=npnts(ix-1)
  160 continue
   ix=ix+1
   if (ix.gt.nmtx) go to 170
   if (mtx(ix).eq.mth) go to 105
   if (mtx(ix).gt.mth) go to 100
   go to 160
  170 continue
   return
   end subroutine rdunf3

   subroutine unresl(ee,tt,sig0,nsig0,sigbkg,sigu,arry)
   !--------------------------------------------------------------------
   ! Computes effective cross-sections in the
   ! unresolved region using the etox method.
   !--------------------------------------------------------------------
   use physics ! provides pi,bk,amassn,amu,hbar,ev
   use util    ! provides error
   use endf    ! provides terpa
   ! externals
   integer::nsig0
   real(kr)::ee,tt,sig0(*),sigbkg(4),sigu(5,*),arry(*)
   ! internals
   integer,parameter::mxks=100
   integer::ntempu,ispot,is,ks,i,isect,mode,inow,nls,lfw
   integer::ne,iest,itp,l,njs,ll,j,ifst,intl,int,is0,kx
   integer::mu,nu,lu,nqf,nqn,nqx,kf,kn,kl,iro,ip,ir,idx
   integer::kkx,ns,knm1,ksp
   real(kr)::k,cwaven,c2,ac,e,spot,sint,abn,awri,ay,spi
   real(kr)::rat,aw,aa,e2,ab,rho,rhoc,vl,ps,aaa
   real(kr)::amux,gxx,amuf,dx,aj,amun,gnox,ggx,gfx,gj,gnx,d
   real(kr)::sigbt,s0u,sti,beta,xj,xk,ttj,yj,yk,ttt
   real(kr)::fact,term,term1,sigtr,enx
   real(kr)::del(3),yy(3),gg(5)
   real(kr),dimension(:,:,:),allocatable,save::sigf
   real(kr),dimension(:,:,:),allocatable,save::tj,tk,tl,t
   real(kr),dimension(:),allocatable,save::abns,sigm
   integer,parameter::nqp=10
   real(kr),dimension(10,4)::qw=reshape((/&
     1.1120413e-1_kr,2.3546798e-1_kr,2.8440987e-1_kr,&
     2.2419127e-1_kr,0.10967668e0_kr,.030493789e0_kr,&
     0.0042930874e0_kr,2.5827047e-4_kr,4.9031965e-6_kr,&
     1.4079206e-8_kr,0.033773418e0_kr,0.079932171e0_kr,&
     0.12835937e0_kr,0.17652616e0_kr,0.21347043e0_kr,&
     0.21154965e0_kr,0.13365186e0_kr,0.022630659e0_kr,&
     1.6313638e-5_kr,2.745383e-31_kr,3.3376214e-4_kr,&
     0.018506108e0_kr,0.12309946e0_kr,0.29918923e0_kr,&
     0.33431475e0_kr,0.17766657e0_kr,0.042695894e0_kr,&
     4.0760575e-3_kr,1.1766115e-4_kr,5.0989546e-7_kr,&
     1.7623788e-3_kr,0.021517749e0_kr,0.080979849e0_kr,&
     0.18797998e0_kr,0.30156335e0_kr,0.29616091e0_kr,&
     0.10775649e0_kr,2.5171914e-3_kr,8.9630388e-10_kr,&
     0.e0_kr/),(/10,4/))
   real(kr),dimension(10,4)::qp=reshape((/&
     3.0013465e-3_kr,7.8592886e-2_kr,4.3282415e-1_kr,&
     1.3345267e0_kr,3.0481846e0_kr,5.8263198e0_kr,&
     9.9452656e0_kr,1.5782128e1_kr,23.996824e0_kr,&
     36.216208e0_kr,1.3219203e-2_kr,7.2349624e-2_kr,&
     0.19089473e0_kr,0.39528842e0_kr,0.74083443e0_kr,&
     1.3498293e0_kr,2.5297983e0_kr,5.2384894e0_kr,&
     13.821772e0_kr,75.647525e0_kr,1.0004488e-3_kr,&
     0.026197629e0_kr,0.14427472e0_kr,0.44484223e0_kr,&
     1.0160615e0_kr,1.9421066e0_kr,3.3150885e0_kr,&
     5.2607092e0_kr,7.9989414e0_kr,12.072069e0_kr,&
     0.013219203e0_kr,0.072349624e0_kr,0.19089473e0_kr,&
     0.39528842e0_kr,0.74083443e0_kr,1.3498293e0_kr,&
     2.5297983e0_kr,5.2384894e0_kr,13.821772e0_kr,&
     75.647525e0_kr/),(/10,4/))
   real(kr),parameter::rc1=.123e0_kr
   real(kr),parameter::rc2=.08e0_kr
   real(kr),parameter::third=.333333333e0_kr
   real(kr),parameter::one=1.e0_kr
   real(kr),parameter::small=1.e-8_kr

   cwaven=sqrt(2*amassn*amu*ev)*1.e-12_kr/hbar

   !--initialize
   c2=bk
   ac=1
   e=abs(ee)
   ntempu=1
   if (allocated(tempu)) deallocate(tempu)
   allocate(tempu(ntempu))
   tempu(1)=tt
   if (tempu(1).lt.one) tempu(1)=one

   !--make a pass to compute the potential scattering cross section
   !--then a second pass to compute unresolved cross sections
   ispot=0
   spot=0
   sint=0
   if (allocated(sigf)) then
      deallocate(sigf)
      deallocate(sigm)
      deallocate(abns)
      deallocate(tj)
      deallocate(tk)
      deallocate(tl)
      deallocate(t)
   endif
   allocate(sigf(ntempu,nsig0,5))
   allocate(sigm(nsig0))
   allocate(t(4,nsig0,ntempu))
   go to 200
  110 continue
   !--use current value of ks to allocate remaining arrays
   allocate(abns(ks))
   allocate(tj(ks,nsig0,ntempu))
   allocate(tk(ks,nsig0,ntempu))
   allocate(tl(3*ks+1,nsig0,ntempu))
   sigbt=sigbkg(1)+spot+sint
   do is=1,nsig0
      sigm(is)=sigbt+sig0(is)
   enddo
   ispot=1

   !--find sections of resonance parameters which contribute
  200 continue
   ks=0
   do 210 i=1,nsect
   isect=i
   abn=abnt(i)
   if (e.lt.elt(i).or.e.gt.eht(i)) go to 210
   mode=modet(i)
   if (mode.lt.5.or.mode.gt.6) go to 210
   inow=ibaset(isect)

   !--retrieve nuclide information
   if (nro.eq.1) then
      iro=inow
      ip=2
      ir=1
      call terpa(ay,e,enx,idx,arry(iro),ip,ir)
      inow=inow+6+2*nint(arry(iro+4))+2*nint(arry(iro+5))
   else
      ay=arry(inow+1)
   endif
   awri=arry(inow)
   aaa=0
   if (naps.eq.2.and.nro.eq.1) aaa=arry(inow+1)
   spi=arry(inow+2)
   if (mode.eq.6) then
      nls=nint(arry(inow+3))
      inow=inow+4
   else
      lfw=nint(arry(inow+3))
      nls=nint(arry(inow+4))
      ne=nint(arry(inow+5))
      inow=inow+6
      if (lfw.ne.0) then
         ! save starting location for fission width energies
         iest=inow
         inow=inow+ne
      endif
   endif

   !--calculate nuclide dependent parameters
   rat=awri/(awri+1)
   aw=awri*amassn
   if (naps.eq.0) then
      aa=rc1*aw**third+rc2
   else if (naps.eq.1) then
      aa=ay
   else if (naps.eq.2.and.nro.eq.1) then
      aa=aaa
   else
      call error('unresl','illegal naps',' ')
   endif
   e2=sqrt(e)
   k=cwaven*rat*e2
   ab=4*pi/k**2
   rho=k*aa
   rhoc=k*ay
   do itp=1,ntempu
      del(itp)=2*e2*sqrt(tempu(itp)*c2/awri)
   enddo

   !--loop over all sequences(l,j) in this section
   do 240 l=1,nls
   njs=nint(arry(inow))
   ll=nint(arry(inow+1))

   !--calculate penetrability and phase shifts
   call uunfac(ll,rho,rhoc,one,vl,ps)
   inow=inow+2
   do 250 j=1,njs

   !--retrieve parameters for this sequence
   ks=ks+1
   if (ks.gt.mxks) call error('unresl',&
                   'storage exceeded, increase mxks',' ')
   if (ispot.ne.0) abns(ks)=abn
   if (mode.ne.6) then
      amux=0
      gxx=0
      if (lfw.eq.0) then
         amuf=0
         dx=arry(inow)
         aj=arry(inow+1)
         amun=arry(inow+2)
         gnox=arry(inow+3)
         ggx=arry(inow+4)
         gfx=0
         inow=inow+5
      else
         amuf=arry(inow)
         dx=arry(inow+1)
         aj=arry(inow+2)
         amun=arry(inow+3)
         gnox=arry(inow+4)
         ggx=arry(inow+5)
         inow=inow+6
         ifst=inow
         intl=2
         call intrf(e,gfx,iest,ifst,ne,intl,arry)
         inow=inow+ne
      endif
   else
      aj=arry(inow)
      int=nint(arry(inow+1))
      int=2
      ne=nint(arry(inow+2))
      amux=arry(inow+3)
      amun=arry(inow+4)
      amuf=arry(inow+6)
      inow=inow+7
      call intr(e,dx,gxx,gnox,ggx,gfx,inow,ne,int,arry)
      inow=inow+6*ne
   endif
   gj=(2*aj+1)/(4*spi+2)
   gnx=gnox*vl*e2*amun

   !--compute potential scattering and interference correction
   if (ispot.ne.1) then
      if (j.eq.1) spot=spot+abn*ab*(2*ll+1)*sin(ps)**2
      sint=sint-abn*ab*pi*gj*gnx*sin(ps)**2/dx

   !--compute etox statistical averages for this sequence
   else
      do is0=1,nsig0
         do itp=1,ntempu
            tk(ks,is0,itp)=0
            do kx=1,4
               t(kx,is0,itp)=0
            enddo
         enddo
      enddo
      d=dx
      mu=nint(amuf)
      if (gfx.le.small) mu=0
      nu=nint(amun)
      lu=nint(amux)
      if (gxx.lt.small) lu=0
      gg(2)=ggx
      if (mu.eq.0) nqf=1
      if (mu.gt.0) nqf=nqp
      if (nu.eq.0) nqn=1
      if (nu.gt.0) nqn=nqp
      if (lu.eq.0) nqx=1
      if (lu.gt.0) nqx=nqp
      do kf=1,nqf
         if (mu.eq.0) gg(1)=gfx
         if (mu.gt.0) gg(1)=qp(kf,mu)*gfx
         do kn=1,nqn
            if (nu.eq.0) gg(3)=gnx
            if (nu.gt.0) gg(3)=qp(kn,nu)*gnx
            do kl=1,nqx
               if (lu.eq.0) gg(4)=gxx
               if (lu.gt.0) gg(4)=qp(kl,lu)*gxx
               gg(5)=gg(1)+gg(2)+gg(3)+gg(4)
               s0u=ab*gj*gg(3)/gg(5)
               do itp=1,ntempu
                  sti=gg(5)/del(itp)
                  do is0=1,nsig0
                     beta=sigm(is0)/s0u
                     call ajku(beta,sti,xj,xk)
                     if (mu.gt.0) xj=xj*qw(kf,mu)
                     if (mu.gt.0) xk=xk*qw(kf,mu)
                     if (nu.gt.0) xj=xj*qw(kn,nu)
                     if (nu.gt.0) xk=xk*qw(kn,nu)
                     if (lu.gt.0) xj=xj*qw(kl,lu)
                     if (lu.gt.0) xk=xk*qw(kl,lu)
                     do kx=1,4
                        t(kx,is0,itp)=t(kx,is0,itp)+xj*gg(kx)
                     enddo
                     tk(ks,is0,itp)=tk(ks,is0,itp)+xk*gg(5)
                  enddo
               enddo
            enddo
         enddo
      enddo
      do itp=1,ntempu
         do is0=1,nsig0
            tk(ks,is0,itp)=ac*tk(ks,is0,itp)/d
            ttj=0
            do kx=1,4
               ttj=ttj+t(kx,is0,itp)
               kkx=kx+(ks-1)*3
               tl(kkx,is0,itp)=ac*t(kx,is0,itp)/d
            enddo
            tj(ks,is0,itp)=ac*ttj/d
         enddo
      enddo
   endif
  250 continue
  240 continue
  210 continue
   if (ispot.eq.0) go to 110

   !--compute average cross-sections by summing over sequences
   ns=ks
   do itp=1,ntempu
      do is0=1,nsig0
         yy(1)=0
         yy(2)=0
         yy(3)=0
         yj=0
         yk=0
         do ks=1,ns
            xk=0
            xj=0
            do ksp=1,ns
               if (ksp.ne.ks) then
                  xj=xj+tj(ksp,is0,itp)*abns(ksp)
                  xk=xk+tk(ksp,is0,itp)*abns(ksp)
               endif
            enddo
            do kx=1,3
               knm1=kx+(ks-1)*3
               yy(kx)=yy(kx)+tl(knm1,is0,itp)*(one-xj)*abns(ks)
            enddo
            ttt=tj(ks,is0,itp)*(one-xj)*abns(ks)
            yj=yj+ttt
            yk=yk+(tk(ks,is0,itp)&
              -tj(ks,is0,itp))*(one-xk)*abns(ks)+ttt
         enddo
         sigf(itp,is0,4)=0
         do i=1,3
            sigf(itp,is0,i)=sigm(is0)*yy(i)/(one-yj)
            sigf(itp,is0,4)=sigf(itp,is0,4)+sigf(itp,is0,i)
         enddo
         fact=sig0(is0)
         term=(one-yj)/(one-yk)
         term1=fact*(yk-yj)/(1-yk)
         fact=sigbt
         sigtr=fact*term+term1
         sigf(itp,is0,5)=sigtr
      enddo
   enddo

   !--output cross section in desired order.
   do is=1,nsig0
      sigu(1,is)=sigf(1,is,4)+sigbt
      sigu(2,is)=sigf(1,is,3)+sigbkg(2)+spot+sint
      sigu(3,is)=sigf(1,is,1)+sigbkg(3)
      sigu(4,is)=sigf(1,is,2)+sigbkg(4)
      sigu(5,is)=sigf(1,is,5)
   enddo
   if (allocated(sigf)) then
      deallocate(sigf)
      deallocate(sigm)
      deallocate(abns)
      deallocate(tj)
      deallocate(tk)
      deallocate(tl)
      deallocate(t)
   endif
   return
   end subroutine unresl

   subroutine uunfac(l,rho,rhoc,amun,vl,ps)
   !-------------------------------------------------------------------
   ! Calculates the penetrability factor (vl) and phase shift (ps).
   !-------------------------------------------------------------------
   ! externals
   integer::l
   real(kr)::rho,rhoc,amun,vl,ps
   ! locals
   real(kr)::r2,r4

   r2=rho*rho
   if (l.eq.0) then
      vl=amun
      ps=rhoc
   else if (l.eq.1) then
      vl=amun*r2/(1+r2)
      vl=amun*r2/(1+r2)
      ps=rhoc-atan(rhoc)
   else
      r4=r2*r2
      vl=amun*r4/(9+3*r2+r4)
      ps=rhoc-atan(3*rhoc/(3-rhoc*rhoc))
   endif
   return
   end subroutine uunfac

   subroutine intrf(e,gfx,iest,ifst,ne,int,a)
   !-------------------------------------------------------------------
   ! Interpolates fission widths for unresolved representation 1.
   !-------------------------------------------------------------------
   use endf ! provides terp1
   ! externals
   integer::iest,ifst,ne,int
   real(kr)::e,gfx
   real(kr)::a(*)
   ! internals
   integer::i,i1,i2

   do i=2,ne
      i1=iest+i-1
      i2=ifst+i-1
      if (e.ge.a(i1-1).and.e.le.a(i1)) then
        call terp1(a(i1-1),a(i2-1),a(i1),a(i2),e,gfx,int)
      endif
   enddo
   return
   end subroutine intrf

   subroutine intr(e,dx,gxx,gnox,ggx,gfx,inow,ne,int,a)
   !-------------------------------------------------------------------
   ! Interpolates energy dependent unresolved resonance parameters.
   !-------------------------------------------------------------------
   use endf ! provides terp1
   ! externals
   integer::inow,ne,int
   real(kr)::e,dx,gxx,gnox,ggx,gfx
   real(kr)::a(*)
   ! internals
   integer::i1,i,i2
   real(kr),parameter::small=1.e-8_kr

   i1=inow
   do i=2,ne
      i2=i1+6
      if (e.ge.a(i1).and.e.le.a(i2)) then
         call terp1(a(i1),a(i1+1),a(i2),a(i2+1),e,dx,int)
         call terp1(a(i1),a(i1+2),a(i2),a(i2+2),e,gxx,int)
         call terp1(a(i1),a(i1+3),a(i2),a(i2+3),e,gnox,int)
         call terp1(a(i1),a(i1+4),a(i2),a(i2+4),e,ggx,int)
         call terp1(a(i1),a(i1+5),a(i2),a(i2+5),e,gfx,int)
      endif
      if (i.lt.ne) i1=i2
   enddo
   if (gxx.lt.small) gxx=0
   if (gfx.lt.small) gfx=0
   return
   end subroutine intr

   subroutine quikw(ax,y,rew,aimw,ki)
   !-------------------------------------------------------------------
   ! Used to calculate chi and xi line shape functions.
   !-------------------------------------------------------------------
   use physics ! provides pi
   ! externals
   integer::ki
   real(kr)::ax,y,rew,aimw
   ! internals
   integer::ii,jj,i,j,n
   real(kr)::rpi,aki,x,test,p,q,p2,q2,pq,hp,hq,hq2,hp2
   real(kr)::a1,a2,a3,a4,a5,d1,d2

   real(kr),parameter::break1=36.e0_kr
   real(kr),parameter::break2=144.e0_kr
   real(kr),parameter::break3=10000.e0_kr
   real(kr),parameter::c1=.2752551e0_kr
   real(kr),parameter::c2=2.724745e0_kr
   real(kr),parameter::c3=.5124242e0_kr
   real(kr),parameter::c4=.05176536e0_kr
   real(kr),parameter::c5=1.1283792e0_kr
   real(kr),parameter::zero=0

   rpi=sqrt(pi)
   aki=1
   if (ax.lt.zero) aki=-1
   x=abs(ax)
   test=x*x+y*y
   if (test.lt.break1) then
      ii=int(x*10)
      jj=int(y*10)
      i=ii+2
      j=jj+2
      n=j-1
      p=10*x-ii
      q=10*y-jj
      p2=p*p
      q2=q*q
      pq=p*q
      hp=p/2
      hq=q/2
      hq2=q2/2
      hp2=p2/2
      a1=hq2-hq
      a2=hp2-hp
      a3=1+pq-p2-q2
      a4=hp2-pq+hp
      a5=hq2-pq+hq
      rew=a1*tr(i,n)+a2*tr(i-1,j)+a3*tr(i,j)+a4*tr(i+1,j)&
        +a5*tr(i,j+1)+pq*tr(i+1,j+1)
      if (ki.gt.0) then
         aimw=a1*ti(i,n)+a2*ti(i-1,j)+a3*ti(i,j)+a4*ti(i+1,j)&
           +a5*ti(i,j+1)+pq*ti(i+1,j+1)
         aimw=aimw*aki
      endif
   else if (test.lt.break2) then
      a1=x**2-y**2
      a2=2*x*y
      a3=a2**2
      a4=a1-c1
      a5=a1-c2
      d1=c3/(a4**2+a3)
      d2=c4/(a5**2+a3)
      rew=d1*(a2*x-a4*y)+d2*(a2*x-a5*y)
      if (ki.gt.0) then
         aimw=d1*(a4*x+a2*y)+d2*(a5*x+a2*y)
         aimw=aimw*aki
      endif
   else if (test.lt.break3) then
      a1=(x**2-y**2)*2
      a2=4*x*y
      a4=a1-1
      d1=c5/(a4**2+a2**2)
      rew=d1*(a2*x-a4*y)
      if (ki.gt.0) then
         aimw=d1*(a4*x+a2*y)
         aimw=aimw*aki
      endif
   else
      a1=1/(rpi*test)
      rew=y*a1
      if (ki.gt.0) aimw=x*a1*aki
   endif
   return
   end subroutine quikw

   subroutine ajku(b,st,aj,ak)
   !-------------------------------------------------------------------
   ! Calculates j and k integrals.
   !-------------------------------------------------------------------
   use physics ! provides pi
   ! externals
   real(kr)::b,st,aj,ak
   ! internals
   integer::n,ki
   real(kr)::pi2,sqpi,term,ep,y,z,yk,zk,b1,a,a2,c,remj,remk
   real(kr)::x1,z1,x2,z2,x3,z3,x4,z4
   real(kr)::y1,ax,rew,aimw
   real(kr),dimension(8)::xg=(/.095012510e0_kr,.28160355e0_kr,&
     .45801678e0_kr,.61787624e0_kr,.75540441e0_kr,.8656312e0_kr,&
     .94457502e0_kr,.98940093e0_kr/)
   real(kr),dimension(8)::wg=(/.18945061e0_kr,.18260342e0_kr,&
     .16915652e0_kr,.14959600e0_kr,.12462897e0_kr,.09515851e0_kr,&
     .06225352e0_kr,.02715246e0_kr/)
   real(kr),parameter::eps=.0002e0_kr
   real(kr),parameter::zero=0

   pi2=pi/2
   sqpi=sqrt(pi)
   if ((-st/2).lt.log(tiny(1._kr))) then
      term=0
   else
      term=exp(-st/2)
   endif
   ep=(1-term)/eps-b
   if (ep.le.zero) then
      aj=pi2/b
      ak=2*aj
   else
      y=0
      z=0
      yk=0
      zk=0
      y1=st/2
      ki=0
      b1=b/(sqpi*y1)
      a=sqrt((1+b)/b)
      a2=a*a
      c=200/(a*st)
      remj=(pi2-atan(c))/(b*a)
      remk=((1+a2)*remj-(1-a2)*c/((c*c+1)*b*a))/(2*a2)
      do n=1,8
         ax=5*xg(n)+5
         call quikw(ax,y1,rew,aimw,ki)
         x1=1/(1+b1/rew)
         z1=x1*x1/rew
         ax=-5*xg(n)+5
         call quikw(ax,y1,rew,aimw,ki)
         x2=1/(1+b1/rew)
         z2=x2*x2/rew
         ax=45*xg(n)+55
         call quikw(ax,y1,rew,aimw,ki)
         x3=1/(1+b1/rew)
         z3=x3*x3/rew
         ax=-45*xg(n)+55
         call quikw(ax,y1,rew,aimw,ki)
         x4=1/(1+b1/rew)
         z4=x4*x4/rew
         y=y+wg(n)*(x3+x4)
         yk=yk+wg(n)*(z3+z4)
         z=z+wg(n)*(x1+x2)
         zk=zk+wg(n)*(z1+z2)
      enddo
      aj=5*(z+9*y)/y1+remj
      ak=aj+5*b1*(zk+9*yk)/y1+remk
   endif
   return
   end subroutine ajku

   subroutine uwtab
   !-------------------------------------------------------------------
   ! Generate table of complex probability integral w.
   !-------------------------------------------------------------------
   ! internals
   integer::nx,ny,nx1,ny1,i,j
   real(kr)::x0,y0,dx,dy,xi,yj,rwt,aimt
   real(kr)::x(62),y(62)
   real(kr),parameter::tenth=0.1e0_kr

   x0=-tenth
   y0=-tenth
   dx=tenth
   dy=tenth
   nx=62
   ny=62
   y(1)=y0
   nx1=nx-1
   ny1=ny-1
   do i=1,ny1
      y(i+1)=i*dy+y0
   enddo
   x(1)=x0
   do i=1,nx1
      x(i+1)=i*dx+x0
   enddo
   do i=1,nx
      do j=1,ny
         xi=x(i)
         yj=y(j)
         call uw(xi,yj,rwt,aimt)
         tr(i,j)=rwt
         ti(i,j)=aimt
      enddo
   enddo
   return
   end subroutine uwtab

   subroutine uw(rez,aim1,rew,aimw)
   !-------------------------------------------------------------------
   ! Compute complex probability integral w.
   !-------------------------------------------------------------------
   use physics ! provides pi
   ! externals
   real(kr)::rez,aim1,rew,aimw
   ! internals
   integer::kw,k,jsig
   real(kr)::rpi,aimz,abrez,aimf
   real(kr)::r2,ai2,rv,ak,el,h,b,a,tempm,temel,g,c,d
   real(kr)::am,aak,ajtemp,temp4,ajp,pr,pim,amagn
   real(kr)::temp1,temp2,aj,ajsig,sigp,expon,expc,exps
   real(kr)::sig2p,aj4sig,aj4sm1,temp3,tt4,temp7,ref
   real(kr)::tempc,tempd
   real(kr),parameter::c1=1.1283792e0_kr
   real(kr),parameter::c2=1.5e0_kr
   real(kr),parameter::up=1.e15_kr
   real(kr),parameter::dn=1.e-15_kr
   real(kr),parameter::brk1=1.25e0_kr
   real(kr),parameter::brk2=5.0e0_kr
   real(kr),parameter::brk3=1.863636e0_kr
   real(kr),parameter::brk4=4.1e0_kr
   real(kr),parameter::brk5=1.71e0_kr
   real(kr),parameter::brk6=2.89e0_kr
   real(kr),parameter::brk7=1.18e0_kr
   real(kr),parameter::brk8=5.76e0_kr
   real(kr),parameter::brk9=1.5e0_kr
   real(kr),parameter::eps=1.e-7_kr
   real(kr),parameter::zero=0

   rpi=sqrt(pi)
   rew=0
   aimw=0
   aimz=abs(aim1)
   abrez=abs(rez)
   if (abrez+aimz.ne.zero) go to 20
   rew=1
   aimw=0
   return

  20 continue
   r2=rez*rez
   ai2=aimz*aimz
   if (abrez+brk1*aimz-brk2.gt.zero) go to 350
   if (abrez+brk3*aimz-brk4.gt.zero) go to 350
   if (r2+brk5*ai2-brk6.lt.zero) go to 340
   if (r2+brk7*ai2-brk8.ge.zero) go to 340
   if (aimz-brk9.ge.zero) go to 350

  340 continue
   kw=2
   aimz=aim1
   go to 420
  350 continue
   kw=1
   if (aim1.ge.zero) go to 370
   kw=2
   aimz=aim1
   go to 420
   ! wa is obtained from asymtotic series
  370 continue
   rv=2*(r2-ai2)
   ak=4*rez*aimz
   el=ak
   h=0
   b=0
   a=0
   tempm=0
   temel=0
   g=1
   c=-c1*aimz
   d=c1*rez
   am=rv-1
   aak=1
   k=0
  380 continue
   ajtemp=2*aak
   temp4=(1-ajtemp)*ajtemp
   ajp=rv-(4*aak+1)
   go to 480
  390 continue
   aak=aak+1
   k=k+1
   pr=rew
   pim=aimw
   amagn=tempm**2+temel**2
   rew=(tempc*tempm+tempd*temel)/amagn
   aimw=(tempm*tempd-temel*tempc)/amagn
   if (abs(rew-pr).ge.eps) go to 380
   if (abs(aimw-pim).ge.eps) go to 380
   return
   ! wt is obtained from taylor series
  420 continue
   temp1=r2+ai2
   temp2=2*temp1*temp1
   aj=-(r2-ai2)/temp2
   ak=2*rez*aimz/temp2
   c=0
   b=0
   ajsig=0
   d=0
   jsig=0
   g=0
   h=0
   el=0
   a=1
   am=1
   sigp=c2
   expon=exp(temp2*aj)
   expc=expon*cos(temp2*ak)
   exps=-expon*sin(temp2*ak)
   sig2p=2*sigp
  430 continue
   aj4sig=4*ajsig
   aj4sm1=aj4sig-1
   temp3=1/(aj4sm1*(aj4sig+3))
   tt4=sig2p*(2*ajsig-1)
   temp4=tt4/(aj4sm1*(aj4sig+1)*(aj4sig-3)*aj4sm1)
   ajp=aj+temp3
   go to 480
  440 continue
   ajsig=ajsig+1
   jsig=jsig+1
   temp7=rpi*(am*am+el*el)
   ref=(aimz*(c*am+d*el)-rez*(am*d-c*el))/temp7/temp1
   aimf=(aimz*(am*d-c*el)+rez*(c*am+d*el))/temp7/temp1
   pr=rew
   pim=aimw
   rew=expc-ref
   aimw=exps-aimf
   if (abs(rew-pr).ge.eps) go to 470
   if (abs(aimw-pim).ge.eps) go to 470
   return
  470 continue
   sig2p=2*ajsig
   go to 430
  480 continue
   tempc=ajp*c+temp4*a-ak*d
   tempd=ajp*d+temp4*b+ak*c
   temel=ajp*el+temp4*h+ak*am
   tempm=ajp*am+temp4*g-ak*el
   a=c
   b=d
   g=am
   h=el
   c=tempc
   d=tempd
   am=tempm
   el=temel
   if ((abs(tempm)+abs(temel)-up).lt.zero) go to 500
   c=dn*c
   d=dn*d
   am=dn*am
   el=dn*el
   tempc=dn*tempc
   tempd=dn*tempd
   tempm=dn*tempm
   temel=dn*temel
   go to 520
  500 continue
   if ((abs(tempm)+abs(temel)-dn).gt.zero) go to 520
   c=up*c
   d=up*d
   am=up*am
   el=up*el
   tempc=up*tempc
   tempd=up*tempd
   tempm=up*tempm
   temel=up*temel
  520 continue
   go to (390,440,530),kw
  530 continue
   return
   end subroutine uw

end module unresm

